// Copyright 2010-2017 Google
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/integral_types.h"
#include "ortools/base/logging.h"
#include "ortools/base/stringprintf.h"
#include "ortools/base/join.h"
#include "ortools/base/map_util.h"
#include "ortools/base/stl_util.h"
#include "ortools/base/mathutil.h"
#include "ortools/constraint_solver/arithmetic.h"
#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"
#include "ortools/util/bitset.h"
#include "ortools/util/saturated_arithmetic.h"
#include "ortools/util/string_array.h"

DEFINE_bool(cp_disable_expression_optimization, false,
            "Disable special optimization when creating expressions.");
DEFINE_bool(cp_share_int_consts, true, "Share IntConst's with the same value.");

#if defined(_MSC_VER)
#pragma warning(disable : 4351 4355)
#endif

namespace operations_research {

namespace arithmetic {

#define COND_REV_ALLOC(rev, alloc) rev ? solver()->RevAlloc(alloc) : alloc;

// ----- Utilities for product expression -----

// Propagates set_min on left * right, left and right >= 0.
void SetPosPosMinExpr(IntExpr* const left, IntExpr* const right, int64 m) {
  DCHECK_GE(left->Min(), 0);
  DCHECK_GE(right->Min(), 0);
  const int64 lmax = left->Max();
  const int64 rmax = right->Max();
  if (m > CapProd(lmax, rmax)) {
    left->solver()->Fail();
  }
  if (m > CapProd(left->Min(), right->Min())) {
    // Ok for m == 0 due to left and right being positive
    if (0 != rmax) {
      left->SetMin(PosIntDivUp(m, rmax));
    }
    if (0 != lmax) {
      right->SetMin(PosIntDivUp(m, lmax));
    }
  }
}

// Propagates set_max on left * right, left and right >= 0.
void SetPosPosMaxExpr(IntExpr* const left, IntExpr* const right, int64 m) {
  DCHECK_GE(left->Min(), 0);
  DCHECK_GE(right->Min(), 0);
  const int64 lmin = left->Min();
  const int64 rmin = right->Min();
  if (m < CapProd(lmin, rmin)) {
    left->solver()->Fail();
  }
  if (m < CapProd(left->Max(), right->Max())) {
    if (0 != lmin) {
      right->SetMax(PosIntDivDown(m, lmin));
    }
    if (0 != rmin) {
      left->SetMax(PosIntDivDown(m, rmin));
    }
    // else do nothing: 0 is supporting any value from other expr.
  }
}

// Propagates set_min on left * right, left >= 0, right across 0.
void SetPosGenMinExpr(IntExpr* const left, IntExpr* const right, int64 m) {
  DCHECK_GE(left->Min(), 0);
  DCHECK_GT(right->Max(), 0);
  DCHECK_LT(right->Min(), 0);
  const int64 lmax = left->Max();
  const int64 rmax = right->Max();
  if (m > CapProd(lmax, rmax)) {
    left->solver()->Fail();
  }
  if (left->Max() == 0) {  // left is bound to 0, product is bound to 0.
    DCHECK_EQ(0, left->Min());
    DCHECK_LE(m, 0);
  } else {
    if (m > 0) {  // We deduce right > 0.
      left->SetMin(PosIntDivUp(m, rmax));
      right->SetMin(PosIntDivUp(m, lmax));
    } else if (m == 0) {
      const int64 lmin = left->Min();
      if (lmin > 0) {
        right->SetMin(0);
      }
    } else {  // m < 0
      const int64 lmin = left->Min();
      if (0 != lmin) {  // We cannot deduce anything if 0 is in the domain.
        right->SetMin(-PosIntDivDown(-m, lmin));
      }
    }
  }
}

// Propagates set_min on left * right, left and right across 0.
void SetGenGenMinExpr(IntExpr* const left, IntExpr* const right, int64 m) {
  DCHECK_LT(left->Min(), 0);
  DCHECK_GT(left->Max(), 0);
  DCHECK_GT(right->Max(), 0);
  DCHECK_LT(right->Min(), 0);
  const int64 lmin = left->Min();
  const int64 lmax = left->Max();
  const int64 rmin = right->Min();
  const int64 rmax = right->Max();
  if (m > std::max(CapProd(lmin, rmin), CapProd(lmax, rmax))) {
    left->solver()->Fail();
  }
  if (m > lmin * rmin) {  // Must be positive section * positive section.
    left->SetMin(PosIntDivUp(m, rmax));
    right->SetMin(PosIntDivUp(m, lmax));
  } else if (m > CapProd(lmax, rmax)) {  // Negative section * negative section.
    left->SetMax(-PosIntDivUp(m, -rmin));
    right->SetMax(-PosIntDivUp(m, -lmin));
  }
}

void TimesSetMin(IntExpr* const left, IntExpr* const right,
                 IntExpr* const minus_left, IntExpr* const minus_right,
                 int64 m) {
  if (left->Min() >= 0) {
    if (right->Min() >= 0) {
      SetPosPosMinExpr(left, right, m);
    } else if (right->Max() <= 0) {
      SetPosPosMaxExpr(left, minus_right, -m);
    } else {  // right->Min() < 0 && right->Max() > 0
      SetPosGenMinExpr(left, right, m);
    }
  } else if (left->Max() <= 0) {
    if (right->Min() >= 0) {
      SetPosPosMaxExpr(right, minus_left, -m);
    } else if (right->Max() <= 0) {
      SetPosPosMinExpr(minus_left, minus_right, m);
    } else {  // right->Min() < 0 && right->Max() > 0
      SetPosGenMinExpr(minus_left, minus_right, m);
    }
  } else if (right->Min() >= 0) {  // left->Min() < 0 && left->Max() > 0
    SetPosGenMinExpr(right, left, m);
  } else if (right->Max() <= 0) {  // left->Min() < 0 && left->Max() > 0
    SetPosGenMinExpr(minus_right, minus_left, m);
  } else {  // left->Min() < 0 && left->Max() > 0 &&
            // right->Min() < 0 && right->Max() > 0
    SetGenGenMinExpr(left, right, m);
  }
}

class TimesIntExpr : public BaseIntExpr {
 public:
  TimesIntExpr(Solver* const s, IntExpr* const l, IntExpr* const r)
      : BaseIntExpr(s),
        left_(l),
        right_(r),
        minus_left_(s->MakeOpposite(left_)),
        minus_right_(s->MakeOpposite(right_)) {}
  ~TimesIntExpr() override {}
  int64 Min() const override {
    const int64 lmin = left_->Min();
    const int64 lmax = left_->Max();
    const int64 rmin = right_->Min();
    const int64 rmax = right_->Max();
    return std::min(std::min(CapProd(lmin, rmin), CapProd(lmax, rmax)),
                    std::min(CapProd(lmax, rmin), CapProd(lmin, rmax)));
  }
  void SetMin(int64 m) override;
  int64 Max() const override {
    const int64 lmin = left_->Min();
    const int64 lmax = left_->Max();
    const int64 rmin = right_->Min();
    const int64 rmax = right_->Max();
    return std::max(std::max(CapProd(lmin, rmin), CapProd(lmax, rmax)),
                    std::max(CapProd(lmax, rmin), CapProd(lmin, rmax)));
  }
  void SetMax(int64 m) override;
  bool Bound() const override;
  std::string name() const override {
    return StringPrintf("(%s * %s)", left_->name().c_str(),
                        right_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s * %s)", left_->DebugString().c_str(),
                        right_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    left_->WhenRange(d);
    right_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kProduct, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, left_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            right_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kProduct, this);
  }

 private:
  IntExpr* const left_;
  IntExpr* const right_;
  IntExpr* const minus_left_;
  IntExpr* const minus_right_;
};

void TimesIntExpr::SetMin(int64 m) {
  if (m != kint64min) {
    TimesSetMin(left_, right_, minus_left_, minus_right_, m);
  }
}

void TimesIntExpr::SetMax(int64 m) {
  if (m != kint64max) {
    TimesSetMin(left_, minus_right_, minus_left_, right_, -m);
  }
}

bool TimesIntExpr::Bound() const {
  const bool left_bound = left_->Bound();
  const bool right_bound = right_->Bound();
  return ((left_bound && left_->Max() == 0) ||
          (right_bound && right_->Max() == 0) || (left_bound && right_bound));
}

// ----- TimesPosIntExpr -----

class TimesPosIntExpr : public BaseIntExpr {
 public:
  TimesPosIntExpr(Solver* const s, IntExpr* const l, IntExpr* const r)
      : BaseIntExpr(s), left_(l), right_(r) {}
  ~TimesPosIntExpr() override {}
  int64 Min() const override { return (left_->Min() * right_->Min()); }
  void SetMin(int64 m) override;
  int64 Max() const override { return (left_->Max() * right_->Max()); }
  void SetMax(int64 m) override;
  bool Bound() const override;
  std::string name() const override {
    return StringPrintf("(%s * %s)", left_->name().c_str(),
                        right_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s * %s)", left_->DebugString().c_str(),
                        right_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    left_->WhenRange(d);
    right_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kProduct, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, left_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            right_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kProduct, this);
  }

 private:
  IntExpr* const left_;
  IntExpr* const right_;
};

void TimesPosIntExpr::SetMin(int64 m) { SetPosPosMinExpr(left_, right_, m); }

void TimesPosIntExpr::SetMax(int64 m) { SetPosPosMaxExpr(left_, right_, m); }

bool TimesPosIntExpr::Bound() const {
  return (left_->Max() == 0 || right_->Max() == 0 ||
          (left_->Bound() && right_->Bound()));
}

// ----- SafeTimesPosIntExpr -----

class SafeTimesPosIntExpr : public BaseIntExpr {
 public:
  SafeTimesPosIntExpr(Solver* const s, IntExpr* const l, IntExpr* const r)
      : BaseIntExpr(s), left_(l), right_(r) {}
  ~SafeTimesPosIntExpr() override {}
  int64 Min() const override { return CapProd(left_->Min(), right_->Min()); }
  void SetMin(int64 m) override {
    if (m != kint64min) {
      SetPosPosMinExpr(left_, right_, m);
    }
  }
  int64 Max() const override { return CapProd(left_->Max(), right_->Max()); }
  void SetMax(int64 m) override {
    if (m != kint64max) {
      SetPosPosMaxExpr(left_, right_, m);
    }
  }
  bool Bound() const override {
    return (left_->Max() == 0 || right_->Max() == 0 ||
            (left_->Bound() && right_->Bound()));
  }
  std::string name() const override {
    return StringPrintf("(%s * %s)", left_->name().c_str(),
                        right_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s * %s)", left_->DebugString().c_str(),
                        right_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    left_->WhenRange(d);
    right_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kProduct, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, left_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            right_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kProduct, this);
  }

 private:
  IntExpr* const left_;
  IntExpr* const right_;
};

// ----- TimesBooleanPosIntExpr -----

class TimesBooleanPosIntExpr : public BaseIntExpr {
 public:
  TimesBooleanPosIntExpr(Solver* const s, BooleanVar* const b, IntExpr* const e)
      : BaseIntExpr(s), boolvar_(b), expr_(e) {}
  ~TimesBooleanPosIntExpr() override {}
  int64 Min() const override {
    return (boolvar_->RawValue() == 1 ? expr_->Min() : 0);
  }
  void SetMin(int64 m) override;
  int64 Max() const override {
    return (boolvar_->RawValue() == 0 ? 0 : expr_->Max());
  }
  void SetMax(int64 m) override;
  void Range(int64* mi, int64* ma) override;
  void SetRange(int64 mi, int64 ma) override;
  bool Bound() const override;
  std::string name() const override {
    return StringPrintf("(%s * %s)", boolvar_->name().c_str(),
                        expr_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s * %s)", boolvar_->DebugString().c_str(),
                        expr_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    boolvar_->WhenRange(d);
    expr_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kProduct, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument,
                                            boolvar_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            expr_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kProduct, this);
  }

 private:
  BooleanVar* const boolvar_;
  IntExpr* const expr_;
};

void TimesBooleanPosIntExpr::SetMin(int64 m) {
  if (m > 0) {
    boolvar_->SetValue(1);
    expr_->SetMin(m);
  }
}

void TimesBooleanPosIntExpr::SetMax(int64 m) {
  if (m < 0) {
    solver()->Fail();
  }
  if (m < expr_->Min()) {
    boolvar_->SetValue(0);
  }
  if (boolvar_->RawValue() == 1) {
    expr_->SetMax(m);
  }
}

void TimesBooleanPosIntExpr::Range(int64* mi, int64* ma) {
  const int value = boolvar_->RawValue();
  if (value == 0) {
    *mi = 0;
    *ma = 0;
  } else if (value == 1) {
    expr_->Range(mi, ma);
  } else {
    *mi = 0;
    *ma = expr_->Max();
  }
}

void TimesBooleanPosIntExpr::SetRange(int64 mi, int64 ma) {
  if (ma < 0 || mi > ma) {
    solver()->Fail();
  }
  if (mi > 0) {
    boolvar_->SetValue(1);
    expr_->SetMin(mi);
  }
  if (ma < expr_->Min()) {
    boolvar_->SetValue(0);
  }
  if (boolvar_->RawValue() == 1) {
    expr_->SetMax(ma);
  }
}

bool TimesBooleanPosIntExpr::Bound() const {
  return (boolvar_->RawValue() == 0 || expr_->Max() == 0 ||
          (boolvar_->RawValue() != BooleanVar::kUnboundBooleanVarValue &&
           expr_->Bound()));
}

// ----- TimesBooleanIntExpr -----

class TimesBooleanIntExpr : public BaseIntExpr {
 public:
  TimesBooleanIntExpr(Solver* const s, BooleanVar* const b, IntExpr* const e)
      : BaseIntExpr(s), boolvar_(b), expr_(e) {}
  ~TimesBooleanIntExpr() override {}
  int64 Min() const override {
    switch (boolvar_->RawValue()) {
      case 0: {
        return 0LL;
      }
      case 1: {
        return expr_->Min();
      }
      default: {
        DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
        return std::min(0LL, expr_->Min());
      }
    }
  }
  void SetMin(int64 m) override;
  int64 Max() const override {
    switch (boolvar_->RawValue()) {
      case 0: {
        return 0LL;
      }
      case 1: {
        return expr_->Max();
      }
      default: {
        DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
        return std::max(0LL, expr_->Max());
      }
    }
  }
  void SetMax(int64 m) override;
  void Range(int64* mi, int64* ma) override;
  void SetRange(int64 mi, int64 ma) override;
  bool Bound() const override;
  std::string name() const override {
    return StringPrintf("(%s * %s)", boolvar_->name().c_str(),
                        expr_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s * %s)", boolvar_->DebugString().c_str(),
                        expr_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    boolvar_->WhenRange(d);
    expr_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kProduct, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument,
                                            boolvar_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            expr_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kProduct, this);
  }

 private:
  BooleanVar* const boolvar_;
  IntExpr* const expr_;
};

void TimesBooleanIntExpr::SetMin(int64 m) {
  switch (boolvar_->RawValue()) {
    case 0: {
      if (m > 0) {
        solver()->Fail();
      }
      break;
    }
    case 1: {
      expr_->SetMin(m);
      break;
    }
    default: {
      DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
      if (m > 0) {  // 0 is no longer possible for boolvar because min > 0.
        boolvar_->SetValue(1);
        expr_->SetMin(m);
      } else if (m <= 0 && expr_->Max() < m) {
        boolvar_->SetValue(0);
      }
    }
  }
}

void TimesBooleanIntExpr::SetMax(int64 m) {
  switch (boolvar_->RawValue()) {
    case 0: {
      if (m < 0) {
        solver()->Fail();
      }
      break;
    }
    case 1: {
      expr_->SetMax(m);
      break;
    }
    default: {
      DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
      if (m < 0) {  // 0 is no longer possible for boolvar because max < 0.
        boolvar_->SetValue(1);
        expr_->SetMax(m);
      } else if (m >= 0 && expr_->Min() > m) {
        boolvar_->SetValue(0);
      }
    }
  }
}

void TimesBooleanIntExpr::Range(int64* mi, int64* ma) {
  switch (boolvar_->RawValue()) {
    case 0: {
      *mi = 0;
      *ma = 0;
      break;
    }
    case 1: {
      *mi = expr_->Min();
      *ma = expr_->Max();
      break;
    }
    default: {
      DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
      *mi = std::min(0LL, expr_->Min());
      *ma = std::max(0LL, expr_->Max());
      break;
    }
  }
}

void TimesBooleanIntExpr::SetRange(int64 mi, int64 ma) {
  if (mi > ma) {
    solver()->Fail();
  }
  switch (boolvar_->RawValue()) {
    case 0: {
      if (mi > 0 || ma < 0) {
        solver()->Fail();
      }
      break;
    }
    case 1: {
      expr_->SetRange(mi, ma);
      break;
    }
    default: {
      DCHECK_EQ(BooleanVar::kUnboundBooleanVarValue, boolvar_->RawValue());
      if (mi > 0) {
        boolvar_->SetValue(1);
        expr_->SetMin(mi);
      } else if (mi == 0 && expr_->Max() < 0) {
        boolvar_->SetValue(0);
      }
      if (ma < 0) {
        boolvar_->SetValue(1);
        expr_->SetMax(ma);
      } else if (ma == 0 && expr_->Min() > 0) {
        boolvar_->SetValue(0);
      }
      break;
    }
  }
}

bool TimesBooleanIntExpr::Bound() const {
  return (boolvar_->RawValue() == 0 ||
          (expr_->Bound() &&
           (boolvar_->RawValue() != BooleanVar::kUnboundBooleanVarValue ||
            expr_->Max() == 0)));
}

// ----- DivPosIntCstExpr -----

class DivPosIntCstExpr : public BaseIntExpr {
 public:
  DivPosIntCstExpr(Solver* const s, IntExpr* const e, int64 v)
      : BaseIntExpr(s), expr_(e), value_(v) {
    CHECK_GE(v, 0);
  }
  ~DivPosIntCstExpr() override {}

  int64 Min() const override { return expr_->Min() / value_; }

  void SetMin(int64 m) override {
    if (m > 0) {
      expr_->SetMin(m * value_);
    } else {
      expr_->SetMin((m - 1) * value_ + 1);
    }
  }
  int64 Max() const override { return expr_->Max() / value_; }

  void SetMax(int64 m) override {
    if (m >= 0) {
      expr_->SetMax((m + 1) * value_ - 1);
    } else {
      expr_->SetMax(m * value_);
    }
  }

  std::string name() const override {
    return StringPrintf("(%s div %" GG_LL_FORMAT "d)", expr_->name().c_str(),
                        value_);
  }

  std::string DebugString() const override {
    return StringPrintf("(%s div %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), value_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kDivide, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kDivide, this);
  }

 private:
  IntExpr* const expr_;
  const int64 value_;
};

// DivPosIntExpr

class DivPosIntExpr : public BaseIntExpr {
 public:
  DivPosIntExpr(Solver* const s, IntExpr* const num, IntExpr* const denom)
      : BaseIntExpr(s),
        num_(num),
        denom_(denom),
        opp_num_(s->MakeOpposite(num)) {}

  ~DivPosIntExpr() override {}

  int64 Min() const override {
    return num_->Min() >= 0
               ? num_->Min() / denom_->Max()
               : (denom_->Min() == 0 ? num_->Min()
                                     : num_->Min() / denom_->Min());
  }

  int64 Max() const override {
    return num_->Max() >= 0 ? (denom_->Min() == 0 ? num_->Max()
                                                  : num_->Max() / denom_->Min())
                            : num_->Max() / denom_->Max();
  }

  static void SetPosMin(IntExpr* const num, IntExpr* const denom, int64 m) {
    num->SetMin(m * denom->Min());
    denom->SetMax(num->Max() / m);
  }

  static void SetPosMax(IntExpr* const num, IntExpr* const denom, int64 m) {
    num->SetMax((m + 1) * denom->Max() - 1);
    denom->SetMin(num->Min() / (m + 1) + 1);
  }

  void SetMin(int64 m) override {
    if (m > 0) {
      SetPosMin(num_, denom_, m);
    } else {
      SetPosMax(opp_num_, denom_, -m);
    }
  }

  void SetMax(int64 m) override {
    if (m >= 0) {
      SetPosMax(num_, denom_, m);
    } else {
      SetPosMin(opp_num_, denom_, -m);
    }
  }

  std::string name() const override {
    return StringPrintf("(%s div %s)", num_->name().c_str(),
                        denom_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s div %s)", num_->DebugString().c_str(),
                        denom_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    num_->WhenRange(d);
    denom_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kDivide, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, num_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            denom_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kDivide, this);
  }

 private:
  IntExpr* const num_;
  IntExpr* const denom_;
  IntExpr* const opp_num_;
};

class DivPosPosIntExpr : public BaseIntExpr {
 public:
  DivPosPosIntExpr(Solver* const s, IntExpr* const num, IntExpr* const denom)
      : BaseIntExpr(s), num_(num), denom_(denom) {}

  ~DivPosPosIntExpr() override {}

  int64 Min() const override {
    if (denom_->Max() == 0) {
      solver()->Fail();
    }
    return num_->Min() / denom_->Max();
  }

  int64 Max() const override {
    if (denom_->Min() == 0) {
      return num_->Max();
    } else {
      return num_->Max() / denom_->Min();
    }
  }

  void SetMin(int64 m) override {
    if (m > 0) {
      num_->SetMin(m * denom_->Min());
      denom_->SetMax(num_->Max() / m);
    }
  }

  void SetMax(int64 m) override {
    if (m >= 0) {
      num_->SetMax((m + 1) * denom_->Max() - 1);
      denom_->SetMin(num_->Min() / (m + 1) + 1);
    } else {
      solver()->Fail();
    }
  }

  std::string name() const override {
    return StringPrintf("(%s div %s)", num_->name().c_str(),
                        denom_->name().c_str());
  }

  std::string DebugString() const override {
    return StringPrintf("(%s div %s)", num_->DebugString().c_str(),
                        denom_->DebugString().c_str());
  }

  void WhenRange(Demon* d) override {
    num_->WhenRange(d);
    denom_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kDivide, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, num_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            denom_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kDivide, this);
  }

 private:
  IntExpr* const num_;
  IntExpr* const denom_;
};

// DivIntExpr

class DivIntExpr : public BaseIntExpr {
 public:
  DivIntExpr(Solver* const s, IntExpr* const num, IntExpr* const denom)
      : BaseIntExpr(s),
        num_(num),
        denom_(denom),
        opp_num_(s->MakeOpposite(num)) {}

  ~DivIntExpr() override {}

  int64 Min() const override {
    const int64 num_min = num_->Min();
    const int64 num_max = num_->Max();
    const int64 denom_min = denom_->Min();
    const int64 denom_max = denom_->Max();

    if (denom_min == 0 && denom_max == 0) {
      return kint64max;  // TODO(user): Check this convention.
    }

    if (denom_min >= 0) {  // Denominator strictly positive.
      DCHECK_GT(denom_max, 0);
      const int64 adjusted_denom_min = denom_min == 0 ? 1 : denom_min;
      return num_min >= 0 ? num_min / denom_max : num_min / adjusted_denom_min;
    } else if (denom_max <= 0) {  // Denominator strictly negative.
      DCHECK_LT(denom_min, 0);
      const int64 adjusted_denom_max = denom_max == 0 ? -1 : denom_max;
      return num_max >= 0 ? num_max / adjusted_denom_max : num_max / denom_min;
    } else {  // Denominator across 0.
      return std::min(num_min, -num_max);
    }
  }

  int64 Max() const override {
    const int64 num_min = num_->Min();
    const int64 num_max = num_->Max();
    const int64 denom_min = denom_->Min();
    const int64 denom_max = denom_->Max();

    if (denom_min == 0 && denom_max == 0) {
      return kint64min;  // TODO(user): Check this convention.
    }

    if (denom_min >= 0) {  // Denominator strictly positive.
      DCHECK_GT(denom_max, 0);
      const int64 adjusted_denom_min = denom_min == 0 ? 1 : denom_min;
      return num_max >= 0 ? num_max / adjusted_denom_min : num_max / denom_max;
    } else if (denom_max <= 0) {  // Denominator strictly negative.
      DCHECK_LT(denom_min, 0);
      const int64 adjusted_denom_max = denom_max == 0 ? -1 : denom_max;
      return num_min >= 0 ? num_min / denom_min
                          : -num_min / -adjusted_denom_max;
    } else {  // Denominator across 0.
      return std::max(num_max, -num_min);
    }
  }

  void AdjustDenominator() {
    if (denom_->Min() == 0) {
      denom_->SetMin(1);
    } else if (denom_->Max() == 0) {
      denom_->SetMax(-1);
    }
  }

  // m > 0.
  static void SetPosMin(IntExpr* const num, IntExpr* const denom, int64 m) {
    DCHECK_GT(m, 0);
    const int64 num_min = num->Min();
    const int64 num_max = num->Max();
    const int64 denom_min = denom->Min();
    const int64 denom_max = denom->Max();
    DCHECK_NE(denom_min, 0);
    DCHECK_NE(denom_max, 0);
    if (denom_min > 0) {  // Denominator strictly positive.
      num->SetMin(m * denom_min);
      denom->SetMax(num_max / m);
    } else if (denom_max < 0) {  // Denominator strictly negative.
      num->SetMax(m * denom_max);
      denom->SetMin(num_min / m);
    } else {  // Denominator across 0.
      if (num_min >= 0) {
        num->SetMin(m);
        denom->SetRange(1, num_max / m);
      } else if (num_max <= 0) {
        num->SetMax(-m);
        denom->SetRange(num_min / m, -1);
      } else {
        if (m > -num_min) {  // Denominator is forced positive.
          num->SetMin(m);
          denom->SetRange(1, num_max / m);
        } else if (m > num_max) {  // Denominator is forced negative.
          num->SetMax(-m);
          denom->SetRange(num_min / m, -1);
        } else {
          denom->SetRange(num_min / m, num_max / m);
        }
      }
    }
  }

  // m >= 0.
  static void SetPosMax(IntExpr* const num, IntExpr* const denom, int64 m) {
    DCHECK_GE(m, 0);
    const int64 num_min = num->Min();
    const int64 num_max = num->Max();
    const int64 denom_min = denom->Min();
    const int64 denom_max = denom->Max();
    DCHECK_NE(denom_min, 0);
    DCHECK_NE(denom_max, 0);
    if (denom_min > 0) {  // Denominator strictly positive.
      num->SetMax((m + 1) * denom_max - 1);
      denom->SetMin((num_min / (m + 1)) + 1);
    } else if (denom_max < 0) {
      num->SetMin((m + 1) * denom_min + 1);
      denom->SetMax(num_max / (m + 1) - 1);
    } else if (num_min > (m + 1) * denom_max - 1) {
      denom->SetMax(-1);
    } else if (num_max < (m + 1) * denom_min + 1) {
      denom->SetMin(1);
    }
  }

  void SetMin(int64 m) override {
    AdjustDenominator();
    if (m > 0) {
      SetPosMin(num_, denom_, m);
    } else {
      SetPosMax(opp_num_, denom_, -m);
    }
  }

  void SetMax(int64 m) override {
    AdjustDenominator();
    if (m >= 0) {
      SetPosMax(num_, denom_, m);
    } else {
      SetPosMin(opp_num_, denom_, -m);
    }
  }

  std::string name() const override {
    return StringPrintf("(%s div %s)", num_->name().c_str(),
                        denom_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("(%s div %s)", num_->DebugString().c_str(),
                        denom_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    num_->WhenRange(d);
    denom_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kDivide, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, num_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            denom_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kDivide, this);
  }

 private:
  IntExpr* const num_;
  IntExpr* const denom_;
  IntExpr* const opp_num_;
};

// ----- IntAbs And IntAbsConstraint ------

class IntAbsConstraint : public CastConstraint {
 public:
  IntAbsConstraint(Solver* const s, IntVar* const sub, IntVar* const target)
      : CastConstraint(s, target), sub_(sub) {}

  ~IntAbsConstraint() override {}

  void Post() override {
    Demon* const sub_demon = MakeConstraintDemon0(
        solver(), this, &IntAbsConstraint::PropagateSub, "PropagateSub");
    sub_->WhenRange(sub_demon);
    Demon* const target_demon = MakeConstraintDemon0(
        solver(), this, &IntAbsConstraint::PropagateTarget, "PropagateTarget");
    target_var_->WhenRange(target_demon);
  }

  void InitialPropagate() override {
    PropagateSub();
    PropagateTarget();
  }

  void PropagateSub() {
    const int64 smin = sub_->Min();
    const int64 smax = sub_->Max();
    if (smax <= 0) {
      target_var_->SetRange(-smax, -smin);
    } else if (smin >= 0) {
      target_var_->SetRange(smin, smax);
    } else {
      target_var_->SetRange(0, std::max(-smin, smax));
    }
  }

  void PropagateTarget() {
    const int64 target_max = target_var_->Max();
    sub_->SetRange(-target_max, target_max);
    const int64 target_min = target_var_->Min();
    if (target_min > 0) {
      if (sub_->Min() > -target_min) {
        sub_->SetMin(target_min);
      } else if (sub_->Max() < target_min) {
        sub_->SetMax(-target_min);
      }
    }
  }

  std::string DebugString() const override {
    return StringPrintf("IntAbsConstraint(%s, %s)", sub_->DebugString().c_str(),
                        target_var_->DebugString().c_str());
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitConstraint(ModelVisitor::kAbsEqual, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            sub_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kTargetArgument,
                                            target_var_);
    visitor->EndVisitConstraint(ModelVisitor::kAbsEqual, this);
  }

 private:
  IntVar* const sub_;
};

class IntAbs : public BaseIntExpr {
 public:
  IntAbs(Solver* const s, IntExpr* const e) : BaseIntExpr(s), expr_(e) {}

  ~IntAbs() override {}

  int64 Min() const override {
    int64 emin = 0;
    int64 emax = 0;
    expr_->Range(&emin, &emax);
    if (emin >= 0) {
      return emin;
    }
    if (emax <= 0) {
      return -emax;
    }
    return 0;
  }

  void SetMin(int64 m) override {
    if (m > 0) {
      int64 emin = 0;
      int64 emax = 0;
      expr_->Range(&emin, &emax);
      if (emin > -m) {
        expr_->SetMin(m);
      } else if (emax < m) {
        expr_->SetMax(-m);
      }
    }
  }

  int64 Max() const override {
    int64 emin = 0;
    int64 emax = 0;
    expr_->Range(&emin, &emax);
    return std::max(-emin, emax);
  }

  void SetMax(int64 m) override { expr_->SetRange(-m, m); }

  void SetRange(int64 mi, int64 ma) override {
    expr_->SetRange(-ma, ma);
    if (mi > 0) {
      int64 emin = 0;
      int64 emax = 0;
      expr_->Range(&emin, &emax);
      if (emin > -mi) {
        expr_->SetMin(mi);
      } else if (emax < mi) {
        expr_->SetMax(-mi);
      }
    }
  }

  void Range(int64* mi, int64* ma) override {
    int64 emin = 0;
    int64 emax = 0;
    expr_->Range(&emin, &emax);
    if (emin >= 0) {
      *mi = emin;
      *ma = emax;
    } else if (emax <= 0) {
      *mi = -emax;
      *ma = -emin;
    } else {
      *mi = 0;
      *ma = std::max(-emin, emax);
    }
  }

  bool Bound() const override { return expr_->Bound(); }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  std::string name() const override {
    return StringPrintf("IntAbs(%s)", expr_->name().c_str());
  }

  std::string DebugString() const override {
    return StringPrintf("IntAbs(%s)", expr_->DebugString().c_str());
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kAbs, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kAbs, this);
  }

  IntVar* CastToVar() override {
    int64 min_value = 0;
    int64 max_value = 0;
    Range(&min_value, &max_value);
    Solver* const s = solver();
    const std::string name = StringPrintf("AbsVar(%s)", expr_->name().c_str());
    IntVar* const target = s->MakeIntVar(min_value, max_value, name);
    CastConstraint* const ct =
        s->RevAlloc(new IntAbsConstraint(s, expr_->Var(), target));
    s->AddCastConstraint(ct, target, this);
    return target;
  }

 private:
  IntExpr* const expr_;
};

// ----- Square -----

// TODO(user): shouldn't we compare to kint32max^2 instead of kint64max?
class IntSquare : public BaseIntExpr {
 public:
  IntSquare(Solver* const s, IntExpr* const e) : BaseIntExpr(s), expr_(e) {}
  ~IntSquare() override {}

  int64 Min() const override {
    const int64 emin = expr_->Min();
    if (emin >= 0) {
      return emin >= kint32max ? kint64max : emin * emin;
    }
    const int64 emax = expr_->Max();
    if (emax < 0) {
      return emax <= -kint32max ? kint64max : emax * emax;
    }
    return 0LL;
  }
  void SetMin(int64 m) override {
    if (m <= 0) {
      return;
    }
    // TODO(user): What happens if m is kint64max?
    const int64 emin = expr_->Min();
    const int64 emax = expr_->Max();
    const int64 root = static_cast<int64>(ceil(sqrt(static_cast<double>(m))));
    if (emin >= 0) {
      expr_->SetMin(root);
    } else if (emax <= 0) {
      expr_->SetMax(-root);
    } else if (expr_->IsVar()) {
      reinterpret_cast<IntVar*>(expr_)->RemoveInterval(-root + 1, root - 1);
    }
  }
  int64 Max() const override {
    const int64 emax = expr_->Max();
    const int64 emin = expr_->Min();
    if (emax >= kint32max || emin <= -kint32max) {
      return kint64max;
    }
    return std::max(emin * emin, emax * emax);
  }
  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m == kint64max) {
      return;
    }
    const int64 root = static_cast<int64>(floor(sqrt(static_cast<double>(m))));
    expr_->SetRange(-root, root);
  }
  bool Bound() const override { return expr_->Bound(); }
  void WhenRange(Demon* d) override { expr_->WhenRange(d); }
  std::string name() const override {
    return StringPrintf("IntSquare(%s)", expr_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("IntSquare(%s)", expr_->DebugString().c_str());
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kSquare, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kSquare, this);
  }

  IntExpr* expr() const { return expr_; }

 protected:
  IntExpr* const expr_;
};

class PosIntSquare : public IntSquare {
 public:
  PosIntSquare(Solver* const s, IntExpr* const e) : IntSquare(s, e) {}
  ~PosIntSquare() override {}

  int64 Min() const override {
    const int64 emin = expr_->Min();
    return emin >= kint32max ? kint64max : emin * emin;
  }
  void SetMin(int64 m) override {
    if (m <= 0) {
      return;
    }
    const int64 root = static_cast<int64>(ceil(sqrt(static_cast<double>(m))));
    expr_->SetMin(root);
  }
  int64 Max() const override {
    const int64 emax = expr_->Max();
    return emax >= kint32max ? kint64max : emax * emax;
  }
  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m == kint64max) {
      return;
    }
    const int64 root = static_cast<int64>(floor(sqrt(static_cast<double>(m))));
    expr_->SetMax(root);
  }
};

// ----- EvenPower -----

int64 IntPower(int64 value, int64 power) {
  int64 result = value;
  // TODO(user): Speed that up.
  for (int i = 1; i < power; ++i) {
    result *= value;
  }
  return result;
}

int64 OverflowLimit(int64 power) {
  return static_cast<int64>(
      floor(exp(log(static_cast<double>(kint64max)) / power)));
}

class BasePower : public BaseIntExpr {
 public:
  BasePower(Solver* const s, IntExpr* const e, int64 n)
      : BaseIntExpr(s), expr_(e), pow_(n), limit_(OverflowLimit(n)) {
    CHECK_GT(n, 0);
  }

  ~BasePower() override {}

  bool Bound() const override { return expr_->Bound(); }

  IntExpr* expr() const { return expr_; }

  int64 exponant() const { return pow_; }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  std::string name() const override {
    return StringPrintf("IntPower(%s, %" GG_LL_FORMAT "d)",
                        expr_->name().c_str(), pow_);
  }

  std::string DebugString() const override {
    return StringPrintf("IntPower(%s, %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), pow_);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kPower, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, pow_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kPower, this);
  }

 protected:
  int64 Pown(int64 value) const {
    if (value >= limit_) {
      return kint64max;
    }
    if (value <= -limit_) {
      if (pow_ % 2 == 0) {
        return kint64max;
      } else {
        return kint64min;
      }
    }
    return IntPower(value, pow_);
  }

  int64 SqrnDown(int64 value) const {
    if (value == kint64min) {
      return kint64min;
    }
    if (value == kint64max) {
      return kint64max;
    }
    int64 res = 0;
    const double d_value = static_cast<double>(value);
    if (value >= 0) {
      const double sq = exp(log(d_value) / pow_);
      res = static_cast<int64>(floor(sq));
    } else {
      CHECK_EQ(1, pow_ % 2);
      const double sq = exp(log(-d_value) / pow_);
      res = -static_cast<int64>(ceil(sq));
    }
    const int64 pow_res = Pown(res + 1);
    if (pow_res <= value) {
      return res + 1;
    } else {
      return res;
    }
  }

  int64 SqrnUp(int64 value) const {
    if (value == kint64min) {
      return kint64min;
    }
    if (value == kint64max) {
      return kint64max;
    }
    int64 res = 0;
    const double d_value = static_cast<double>(value);
    if (value >= 0) {
      const double sq = exp(log(d_value) / pow_);
      res = static_cast<int64>(ceil(sq));
    } else {
      CHECK_EQ(1, pow_ % 2);
      const double sq = exp(log(-d_value) / pow_);
      res = -static_cast<int64>(floor(sq));
    }
    const int64 pow_res = Pown(res - 1);
    if (pow_res >= value) {
      return res - 1;
    } else {
      return res;
    }
  }

  IntExpr* const expr_;
  const int64 pow_;
  const int64 limit_;
};

class IntEvenPower : public BasePower {
 public:
  IntEvenPower(Solver* const s, IntExpr* const e, int64 n)
      : BasePower(s, e, n) {
    CHECK_EQ(0, n % 2);
  }

  ~IntEvenPower() override {}

  int64 Min() const override {
    int64 emin = 0;
    int64 emax = 0;
    expr_->Range(&emin, &emax);
    if (emin >= 0) {
      return Pown(emin);
    }
    if (emax < 0) {
      return Pown(emax);
    }
    return 0LL;
  }
  void SetMin(int64 m) override {
    if (m <= 0) {
      return;
    }
    int64 emin = 0;
    int64 emax = 0;
    expr_->Range(&emin, &emax);
    const int64 root = SqrnUp(m);
    if (emin > -root) {
      expr_->SetMin(root);
    } else if (emax < root) {
      expr_->SetMax(-root);
    } else if (expr_->IsVar()) {
      reinterpret_cast<IntVar*>(expr_)->RemoveInterval(-root + 1, root - 1);
    }
  }

  int64 Max() const override {
    return std::max(Pown(expr_->Min()), Pown(expr_->Max()));
  }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m == kint64max) {
      return;
    }
    const int64 root = SqrnDown(m);
    expr_->SetRange(-root, root);
  }
};

class PosIntEvenPower : public BasePower {
 public:
  PosIntEvenPower(Solver* const s, IntExpr* const e, int64 pow)
      : BasePower(s, e, pow) {
    CHECK_EQ(0, pow % 2);
  }

  ~PosIntEvenPower() override {}

  int64 Min() const override { return Pown(expr_->Min()); }

  void SetMin(int64 m) override {
    if (m <= 0) {
      return;
    }
    expr_->SetMin(SqrnUp(m));
  }
  int64 Max() const override { return Pown(expr_->Max()); }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m == kint64max) {
      return;
    }
    expr_->SetMax(SqrnDown(m));
  }
};

class IntOddPower : public BasePower {
 public:
  IntOddPower(Solver* const s, IntExpr* const e, int64 n) : BasePower(s, e, n) {
    CHECK_EQ(1, n % 2);
  }

  ~IntOddPower() override {}

  int64 Min() const override { return Pown(expr_->Min()); }

  void SetMin(int64 m) override { expr_->SetMin(SqrnUp(m)); }

  int64 Max() const override { return Pown(expr_->Max()); }

  void SetMax(int64 m) override { expr_->SetMax(SqrnDown(m)); }
};

// ----- Min(expr, expr) -----

class MinIntExpr : public BaseIntExpr {
 public:
  MinIntExpr(Solver* const s, IntExpr* const l, IntExpr* const r)
      : BaseIntExpr(s), left_(l), right_(r) {}
  ~MinIntExpr() override {}
  int64 Min() const override {
    const int64 lmin = left_->Min();
    const int64 rmin = right_->Min();
    return std::min(lmin, rmin);
  }
  void SetMin(int64 m) override {
    left_->SetMin(m);
    right_->SetMin(m);
  }
  int64 Max() const override {
    const int64 lmax = left_->Max();
    const int64 rmax = right_->Max();
    return std::min(lmax, rmax);
  }
  void SetMax(int64 m) override {
    if (left_->Min() > m) {
      right_->SetMax(m);
    }
    if (right_->Min() > m) {
      left_->SetMax(m);
    }
  }
  std::string name() const override {
    return StringPrintf("MinIntExpr(%s, %s)", left_->name().c_str(),
                        right_->name().c_str());
  }
  std::string DebugString() const override {
    return StringPrintf("MinIntExpr(%s, %s)", left_->DebugString().c_str(),
                        right_->DebugString().c_str());
  }
  void WhenRange(Demon* d) override {
    left_->WhenRange(d);
    right_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kMin, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, left_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            right_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kMin, this);
  }

 private:
  IntExpr* const left_;
  IntExpr* const right_;
};

// ----- Min(expr, constant) -----

class MinCstIntExpr : public BaseIntExpr {
 public:
  MinCstIntExpr(Solver* const s, IntExpr* const e, int64 v)
      : BaseIntExpr(s), expr_(e), value_(v) {}

  ~MinCstIntExpr() override {}

  int64 Min() const override { return std::min(expr_->Min(), value_); }

  void SetMin(int64 m) override {
    if (m > value_) {
      solver()->Fail();
    }
    expr_->SetMin(m);
  }

  int64 Max() const override { return std::min(expr_->Max(), value_); }

  void SetMax(int64 m) override {
    if (value_ > m) {
      expr_->SetMax(m);
    }
  }

  bool Bound() const override {
    return (expr_->Bound() || expr_->Min() >= value_);
  }

  std::string name() const override {
    return StringPrintf("MinCstIntExpr(%s, %" GG_LL_FORMAT "d)",
                        expr_->name().c_str(), value_);
  }

  std::string DebugString() const override {
    return StringPrintf("MinCstIntExpr(%s, %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), value_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kMin, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kMin, this);
  }

 private:
  IntExpr* const expr_;
  const int64 value_;
};

// ----- Max(expr, expr) -----

class MaxIntExpr : public BaseIntExpr {
 public:
  MaxIntExpr(Solver* const s, IntExpr* const l, IntExpr* const r)
      : BaseIntExpr(s), left_(l), right_(r) {}

  ~MaxIntExpr() override {}

  int64 Min() const override { return std::max(left_->Min(), right_->Min()); }

  void SetMin(int64 m) override {
    if (left_->Max() < m) {
      right_->SetMin(m);
    } else {
      if (right_->Max() < m) {
        left_->SetMin(m);
      }
    }
  }

  int64 Max() const override { return std::max(left_->Max(), right_->Max()); }

  void SetMax(int64 m) override {
    left_->SetMax(m);
    right_->SetMax(m);
  }

  std::string name() const override {
    return StringPrintf("MaxIntExpr(%s, %s)", left_->name().c_str(),
                        right_->name().c_str());
  }

  std::string DebugString() const override {
    return StringPrintf("MaxIntExpr(%s, %s)", left_->DebugString().c_str(),
                        right_->DebugString().c_str());
  }

  void WhenRange(Demon* d) override {
    left_->WhenRange(d);
    right_->WhenRange(d);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kMax, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kLeftArgument, left_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kRightArgument,
                                            right_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kMax, this);
  }

 private:
  IntExpr* const left_;
  IntExpr* const right_;
};

// ----- Max(expr, constant) -----

class MaxCstIntExpr : public BaseIntExpr {
 public:
  MaxCstIntExpr(Solver* const s, IntExpr* const e, int64 v)
      : BaseIntExpr(s), expr_(e), value_(v) {}

  ~MaxCstIntExpr() override {}

  int64 Min() const override { return std::max(expr_->Min(), value_); }

  void SetMin(int64 m) override {
    if (value_ < m) {
      expr_->SetMin(m);
    }
  }

  int64 Max() const override { return std::max(expr_->Max(), value_); }

  void SetMax(int64 m) override {
    if (m < value_) {
      solver()->Fail();
    }
    expr_->SetMax(m);
  }

  bool Bound() const override {
    return (expr_->Bound() || expr_->Max() <= value_);
  }

  std::string name() const override {
    return StringPrintf("MaxCstIntExpr(%s, %" GG_LL_FORMAT "d)",
                        expr_->name().c_str(), value_);
  }

  std::string DebugString() const override {
    return StringPrintf("MaxCstIntExpr(%s, %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), value_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kMax, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kMax, this);
  }

 private:
  IntExpr* const expr_;
  const int64 value_;
};

// ----- Convex Piecewise -----

// This class is a very simple convex piecewise linear function.  The
// argument of the function is the expression.  Between early_date and
// late_date, the value of the function is 0.  Before early date, it
// is affine and the cost is early_cost * (early_date - x). After
// late_date, the cost is late_cost * (x - late_date).

class SimpleConvexPiecewiseExpr : public BaseIntExpr {
 public:
  SimpleConvexPiecewiseExpr(Solver* const s, IntExpr* const e, int64 ec,
                            int64 ed, int64 ld, int64 lc)
      : BaseIntExpr(s),
        expr_(e),
        early_cost_(ec),
        early_date_(ec == 0 ? kint64min : ed),
        late_date_(lc == 0 ? kint64max : ld),
        late_cost_(lc) {
    DCHECK_GE(ec, 0LL);
    DCHECK_GE(lc, 0LL);
    DCHECK_GE(ld, ed);

    // If the penalty is 0, we can push the "confort zone or zone
    // of no cost towards infinity.
  }

  ~SimpleConvexPiecewiseExpr() override {}

  int64 Min() const override {
    const int64 vmin = expr_->Min();
    const int64 vmax = expr_->Max();
    if (vmin >= late_date_) {
      return (vmin - late_date_) * late_cost_;
    } else if (vmax <= early_date_) {
      return (early_date_ - vmax) * early_cost_;
    } else {
      return 0LL;
    }
  }

  void SetMin(int64 m) override {
    if (m <= 0) {
      return;
    }
    int64 vmin = 0;
    int64 vmax = 0;
    expr_->Range(&vmin, &vmax);

    const int64 rb =
        (late_cost_ == 0 ? vmax : late_date_ + PosIntDivUp(m, late_cost_) - 1);
    const int64 lb =
        (early_cost_ == 0 ? vmin
                          : early_date_ - PosIntDivUp(m, early_cost_) + 1);

    if (expr_->IsVar()) {
      expr_->Var()->RemoveInterval(lb, rb);
    }
  }

  int64 Max() const override {
    const int64 vmin = expr_->Min();
    const int64 vmax = expr_->Max();
    const int64 mr = vmax > late_date_ ? (vmax - late_date_) * late_cost_ : 0;
    const int64 ml =
        vmin < early_date_ ? (early_date_ - vmin) * early_cost_ : 0;
    return std::max(mr, ml);
  }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (late_cost_ != 0LL) {
      const int64 rb = late_date_ + PosIntDivDown(m, late_cost_);
      if (early_cost_ != 0LL) {
        const int64 lb = early_date_ - PosIntDivDown(m, early_cost_);
        expr_->SetRange(lb, rb);
      } else {
        expr_->SetMax(rb);
      }
    } else {
      if (early_cost_ != 0LL) {
        const int64 lb = early_date_ - PosIntDivDown(m, early_cost_);
        expr_->SetMin(lb);
      }
    }
  }

  std::string name() const override {
    return StringPrintf("ConvexPiecewiseExpr(%s, ec = %" GG_LL_FORMAT
                        "d, ed = %" GG_LL_FORMAT "d, ld = %" GG_LL_FORMAT
                        "d, lc = %" GG_LL_FORMAT "d)",
                        expr_->name().c_str(), early_cost_, early_date_,
                        late_date_, late_cost_);
  }

  std::string DebugString() const override {
    return StringPrintf("ConvexPiecewiseExpr(%s, ec = %" GG_LL_FORMAT
                        "d, ed = %" GG_LL_FORMAT "d, ld = %" GG_LL_FORMAT
                        "d, lc = %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), early_cost_, early_date_,
                        late_date_, late_cost_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kConvexPiecewise, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kEarlyCostArgument,
                                  early_cost_);
    visitor->VisitIntegerArgument(ModelVisitor::kEarlyDateArgument,
                                  early_date_);
    visitor->VisitIntegerArgument(ModelVisitor::kLateCostArgument, late_cost_);
    visitor->VisitIntegerArgument(ModelVisitor::kLateDateArgument, late_date_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kConvexPiecewise, this);
  }

 private:
  IntExpr* const expr_;
  const int64 early_cost_;
  const int64 early_date_;
  const int64 late_date_;
  const int64 late_cost_;
};

// ----- Semi Continuous -----

class SemiContinuousExpr : public BaseIntExpr {
 public:
  SemiContinuousExpr(Solver* const s, IntExpr* const e, int64 fixed_charge,
                     int64 step)
      : BaseIntExpr(s), expr_(e), fixed_charge_(fixed_charge), step_(step) {
    DCHECK_GE(fixed_charge, 0LL);
    DCHECK_GT(step, 0LL);
  }

  ~SemiContinuousExpr() override {}

  int64 Value(int64 x) const {
    if (x <= 0) {
      return 0;
    } else {
      return CapAdd(fixed_charge_, CapProd(x, step_));
    }
  }

  int64 Min() const override { return Value(expr_->Min()); }

  void SetMin(int64 m) override {
    if (m >= CapAdd(fixed_charge_, step_)) {
      const int64 y = PosIntDivUp(CapSub(m, fixed_charge_), step_);
      expr_->SetMin(y);
    } else if (m > 0) {
      expr_->SetMin(1);
    }
  }

  int64 Max() const override { return Value(expr_->Max()); }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m == kint64max) {
      return;
    }
    if (m < CapAdd(fixed_charge_, step_)) {
      expr_->SetMax(0);
    } else {
      const int64 y = PosIntDivDown(CapSub(m, fixed_charge_), step_);
      expr_->SetMax(y);
    }
  }

  std::string name() const override {
    return StringPrintf("SemiContinuous(%s, fixed_charge = %" GG_LL_FORMAT
                        "d, step = %" GG_LL_FORMAT "d)",
                        expr_->name().c_str(), fixed_charge_, step_);
  }

  std::string DebugString() const override {
    return StringPrintf("SemiContinuous(%s, fixed_charge = %" GG_LL_FORMAT
                        "d, step = %" GG_LL_FORMAT "d)",
                        expr_->DebugString().c_str(), fixed_charge_, step_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kFixedChargeArgument,
                                  fixed_charge_);
    visitor->VisitIntegerArgument(ModelVisitor::kStepArgument, step_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
  }

 private:
  IntExpr* const expr_;
  const int64 fixed_charge_;
  const int64 step_;
};

class SemiContinuousStepOneExpr : public BaseIntExpr {
 public:
  SemiContinuousStepOneExpr(Solver* const s, IntExpr* const e,
                            int64 fixed_charge)
      : BaseIntExpr(s), expr_(e), fixed_charge_(fixed_charge) {
    DCHECK_GE(fixed_charge, 0LL);
  }

  ~SemiContinuousStepOneExpr() override {}

  int64 Value(int64 x) const {
    if (x <= 0) {
      return 0;
    } else {
      return fixed_charge_ + x;
    }
  }

  int64 Min() const override { return Value(expr_->Min()); }

  void SetMin(int64 m) override {
    if (m >= fixed_charge_ + 1) {
      expr_->SetMin(m - fixed_charge_);
    } else if (m > 0) {
      expr_->SetMin(1);
    }
  }

  int64 Max() const override { return Value(expr_->Max()); }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m < fixed_charge_ + 1) {
      expr_->SetMax(0);
    } else {
      expr_->SetMax(m - fixed_charge_);
    }
  }

  std::string name() const override {
    return StringPrintf(
        "SemiContinuousStepOne(%s, fixed_charge = %" GG_LL_FORMAT "d)",
        expr_->name().c_str(), fixed_charge_);
  }

  std::string DebugString() const override {
    return StringPrintf(
        "SemiContinuousStepOne(%s, fixed_charge = %" GG_LL_FORMAT "d)",
        expr_->DebugString().c_str(), fixed_charge_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kFixedChargeArgument,
                                  fixed_charge_);
    visitor->VisitIntegerArgument(ModelVisitor::kStepArgument, 1);
    visitor->EndVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
  }

 private:
  IntExpr* const expr_;
  const int64 fixed_charge_;
};

class SemiContinuousStepZeroExpr : public BaseIntExpr {
 public:
  SemiContinuousStepZeroExpr(Solver* const s, IntExpr* const e,
                             int64 fixed_charge)
      : BaseIntExpr(s), expr_(e), fixed_charge_(fixed_charge) {
    DCHECK_GT(fixed_charge, 0LL);
  }

  ~SemiContinuousStepZeroExpr() override {}

  int64 Value(int64 x) const {
    if (x <= 0) {
      return 0;
    } else {
      return fixed_charge_;
    }
  }

  int64 Min() const override { return Value(expr_->Min()); }

  void SetMin(int64 m) override {
    if (m >= fixed_charge_) {
      solver()->Fail();
    } else if (m > 0) {
      expr_->SetMin(1);
    }
  }

  int64 Max() const override { return Value(expr_->Max()); }

  void SetMax(int64 m) override {
    if (m < 0) {
      solver()->Fail();
    }
    if (m < fixed_charge_) {
      expr_->SetMax(0);
    }
  }

  std::string name() const override {
    return StringPrintf(
        "SemiContinuousStepZero(%s, fixed_charge = %" GG_LL_FORMAT "d)",
        expr_->name().c_str(), fixed_charge_);
  }

  std::string DebugString() const override {
    return StringPrintf(
        "SemiContinuousStepZero(%s, fixed_charge = %" GG_LL_FORMAT "d)",
        expr_->DebugString().c_str(), fixed_charge_);
  }

  void WhenRange(Demon* d) override { expr_->WhenRange(d); }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerArgument(ModelVisitor::kFixedChargeArgument,
                                  fixed_charge_);
    visitor->VisitIntegerArgument(ModelVisitor::kStepArgument, 0);
    visitor->EndVisitIntegerExpression(ModelVisitor::kSemiContinuous, this);
  }

 private:
  IntExpr* const expr_;
  const int64 fixed_charge_;
};

// This constraints links an expression and the variable it is casted into
class LinkExprAndVar : public CastConstraint {
 public:
  LinkExprAndVar(Solver* const s, IntExpr* const expr, IntVar* const var)
      : CastConstraint(s, var), expr_(expr) {}

  ~LinkExprAndVar() override {}

  void Post() override {
    Solver* const s = solver();
    Demon* d = s->MakeConstraintInitialPropagateCallback(this);
    expr_->WhenRange(d);
    target_var_->WhenRange(d);
  }

  void InitialPropagate() override {
    expr_->SetRange(target_var_->Min(), target_var_->Max());
    int64 l, u;
    expr_->Range(&l, &u);
    target_var_->SetRange(l, u);
  }

  std::string DebugString() const override {
    return StringPrintf("cast(%s, %s)", expr_->DebugString().c_str(),
                        target_var_->DebugString().c_str());
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitConstraint(ModelVisitor::kLinkExprVar, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expr_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kTargetArgument,
                                            target_var_);
    visitor->EndVisitConstraint(ModelVisitor::kLinkExprVar, this);
  }

 private:
  IntExpr* const expr_;
};

// ----- Conditional Expression -----

class ExprWithEscapeValue : public BaseIntExpr {
 public:
  ExprWithEscapeValue(Solver* const s, IntVar* const c, IntExpr* const e,
                      int64 unperformed_value)
      : BaseIntExpr(s),
        condition_(c),
        expression_(e),
        unperformed_value_(unperformed_value) {}

  ~ExprWithEscapeValue() override {}

  int64 Min() const override {
    if (condition_->Min() == 1) {
      return expression_->Min();
    } else if (condition_->Max() == 1) {
      return std::min(unperformed_value_, expression_->Min());
    } else {
      return unperformed_value_;
    }
  }

  void SetMin(int64 m) override {
    if (m > unperformed_value_) {
      condition_->SetValue(1);
      expression_->SetMin(m);
    } else if (condition_->Min() == 1) {
      expression_->SetMin(m);
    } else if (m > expression_->Max()) {
      condition_->SetValue(0);
    }
  }

  int64 Max() const override {
    if (condition_->Min() == 1) {
      return expression_->Max();
    } else if (condition_->Max() == 1) {
      return std::max(unperformed_value_, expression_->Max());
    } else {
      return unperformed_value_;
    }
  }

  void SetMax(int64 m) override {
    if (m < unperformed_value_) {
      condition_->SetValue(1);
      expression_->SetMax(m);
    } else if (condition_->Min() == 1) {
      expression_->SetMax(m);
    } else if (m < expression_->Min()) {
      condition_->SetValue(0);
    }
  }

  void SetRange(int64 mi, int64 ma) override {
    if (ma < unperformed_value_ || mi > unperformed_value_) {
      condition_->SetValue(1);
      expression_->SetRange(mi, ma);
    } else if (condition_->Min() == 1) {
      expression_->SetRange(mi, ma);
    } else if (ma < expression_->Min() || mi > expression_->Max()) {
      condition_->SetValue(0);
    }
  }

  void SetValue(int64 v) override {
    if (v != unperformed_value_) {
      condition_->SetValue(1);
      expression_->SetValue(v);
    } else if (condition_->Min() == 1) {
      expression_->SetValue(v);
    } else if (v < expression_->Min() || v > expression_->Max()) {
      condition_->SetValue(0);
    }
  }

  bool Bound() const override {
    return condition_->Max() == 0 || expression_->Bound();
  }

  void WhenRange(Demon* d) override {
    expression_->WhenRange(d);
    condition_->WhenBound(d);
  }

  std::string DebugString() const override {
    return StringPrintf("ConditionExpr(%s, %s, %" GG_LL_FORMAT "d)",
                        condition_->DebugString().c_str(),
                        expression_->DebugString().c_str(), unperformed_value_);
  }

  void Accept(ModelVisitor* const visitor) const override {
    visitor->BeginVisitIntegerExpression(ModelVisitor::kConditionalExpr, this);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kVariableArgument,
                                            condition_);
    visitor->VisitIntegerExpressionArgument(ModelVisitor::kExpressionArgument,
                                            expression_);
    visitor->VisitIntegerArgument(ModelVisitor::kValueArgument,
                                  unperformed_value_);
    visitor->EndVisitIntegerExpression(ModelVisitor::kConditionalExpr, this);
  }

 private:
  IntVar* const condition_;
  IntExpr* const expression_;
  const int64 unperformed_value_;
  DISALLOW_COPY_AND_ASSIGN(ExprWithEscapeValue);
};

void ExtractPower(IntExpr** const expr, int64* const exponant) {
  if (dynamic_cast<BasePower*>(*expr) != nullptr) {
    BasePower* const power = dynamic_cast<BasePower*>(*expr);
    *expr = power->expr();
    *exponant = power->exponant();
  }
  if (dynamic_cast<IntSquare*>(*expr) != nullptr) {
    IntSquare* const power = dynamic_cast<IntSquare*>(*expr);
    *expr = power->expr();
    *exponant = 2;
  }
  if ((*expr)->IsVar()) {
    IntVar* const var = (*expr)->Var();
    IntExpr* const sub = var->solver()->CastExpression(var);
    if (sub != nullptr && dynamic_cast<BasePower*>(sub) != nullptr) {
      BasePower* const power = dynamic_cast<BasePower*>(sub);
      *expr = power->expr();
      *exponant = power->exponant();
    }
    if (sub != nullptr && dynamic_cast<IntSquare*>(sub) != nullptr) {
      IntSquare* const power = dynamic_cast<IntSquare*>(sub);
      *expr = power->expr();
      *exponant = 2;
    }
  }
}

#undef COND_REV_ALLOC
}


}  // namespace operations_research
