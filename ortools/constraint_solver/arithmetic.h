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

#ifndef OR_TOOLS_CONSTRAINT_SOLVER_ARITHMETIC_H_
#define OR_TOOLS_CONSTRAINT_SOLVER_ARITHMETIC_H_

#if defined(_MSC_VER)
#pragma warning(disable : 4351 4355)
#endif

#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"

namespace operations_research {

namespace arithmetic {

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

class UnaryIterator : public IntVarIterator {
 public:
  UnaryIterator(const IntVar* const v, bool hole, bool reversible)
      : iterator_(hole ? v->MakeHoleIterator(reversible)
                       : v->MakeDomainIterator(reversible)),
        reversible_(reversible) {}

  ~UnaryIterator() override {
    if (!reversible_) {
      delete iterator_;
    }
  }

  void Init() override { iterator_->Init(); }

  bool Ok() const override { return iterator_->Ok(); }

  void Next() override { iterator_->Next(); }

 protected:
  IntVarIterator* const iterator_;
  const bool reversible_;
};

}

}  // namespace operations_research

#endif  // OR_TOOLS_CONSTRAINT_SOLVER_ARITHMETIC_H_
